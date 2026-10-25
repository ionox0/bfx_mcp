"""Shared HTTP layer: one async client, an on-disk response cache, per-host
pacing, and retries.

All three upstreams are free public services with rate limits (NCBI: 3 req/s
without a key, Ensembl: ~15 req/s, gnomAD: unpublished but strict), and Ensembl
in particular returns transient 503s and slow responses. Caching every
successful response on disk keeps us polite, makes repeated agent calls cheap,
and makes eval runs reproducible.
"""

import asyncio
import hashlib
import json
import logging
import os
import random
import sqlite3
import sys
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import httpx2 as httpx

USER_AGENT = "genomics-mcp/0.1 (+https://modelcontextprotocol.io)"
CACHE_TTL_S = float(os.environ.get("GENOMICS_MCP_CACHE_TTL_S", 7 * 24 * 3600))
CACHE_PATH = Path(os.environ.get("GENOMICS_MCP_CACHE", Path.home() / ".cache" / "genomics-mcp" / "http.sqlite"))

# Minimum seconds between requests to the same host.
_MIN_INTERVAL = {
    "eutils.ncbi.nlm.nih.gov": 0.11 if os.environ.get("NCBI_API_KEY") else 0.35,
    "rest.ensembl.org": 0.07,
    "gnomad.broadinstitute.org": 1.0,
    "api.genome.ucsc.edu": 1.0,
}
_RETRY_STATUS = {429, 500, 502, 503, 504}

# The HTTP library logs every request at INFO; retries are logged below instead.
for _name in ("httpx", "httpx2"):
    logging.getLogger(_name).setLevel(logging.WARNING)


class UpstreamError(Exception):
    """An upstream API failed after retries, or returned a non-retryable error."""


def log(msg: str) -> None:
    # stdout is the MCP channel under stdio, so diagnostics go to stderr.
    print(f"[genomics-mcp] {msg}", file=sys.stderr, flush=True)


class _Cache:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        # timeout: the eval runner and dataset builder may share one cache file.
        self._db = sqlite3.connect(path, check_same_thread=False, timeout=30)
        self._db.execute("CREATE TABLE IF NOT EXISTS http (key TEXT PRIMARY KEY, ts REAL, body TEXT)")

    def get(self, key: str) -> str | None:
        row = self._db.execute("SELECT ts, body FROM http WHERE key = ?", (key,)).fetchone()
        if row and time.time() - row[0] < CACHE_TTL_S:
            return row[1]
        return None

    def put(self, key: str, body: str) -> None:
        self._db.execute("INSERT OR REPLACE INTO http VALUES (?, ?, ?)", (key, time.time(), body))
        self._db.commit()


class Http:
    def __init__(self) -> None:
        self._client = httpx.AsyncClient(
            timeout=httpx.Timeout(60.0, connect=15.0),
            headers={"User-Agent": USER_AGENT},
            follow_redirects=True,
        )
        self._cache = _Cache(CACHE_PATH)
        self._locks: dict[str, asyncio.Lock] = {}
        self._last: dict[str, float] = {}

    async def _pace(self, host: str) -> None:
        lock = self._locks.setdefault(host, asyncio.Lock())
        async with lock:
            wait = self._last.get(host, 0) + _MIN_INTERVAL.get(host, 0.1) - time.monotonic()
            if wait > 0:
                await asyncio.sleep(wait)
            self._last[host] = time.monotonic()

    async def request(
        self,
        method: str,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        json_body: Any = None,
        headers: dict[str, str] | None = None,
        attempts: int = 4,
    ) -> str:
        """Return the response body as text, served from cache when fresh."""
        key = hashlib.sha256(
            json.dumps([method, url, params, json_body], sort_keys=True, default=str).encode()
        ).hexdigest()
        if (hit := self._cache.get(key)) is not None:
            return hit

        host = urlsplit(url).hostname or ""
        last_err = ""
        for attempt in range(attempts):
            await self._pace(host)
            try:
                resp = await self._client.request(method, url, params=params, json=json_body, headers=headers)
            except (httpx.TransportError, httpx.TimeoutException) as e:
                last_err = f"{type(e).__name__}: {e}"
            else:
                if resp.status_code < 400:
                    self._cache.put(key, resp.text)
                    return resp.text
                last_err = f"HTTP {resp.status_code}: {resp.text[:300]}"
                if resp.status_code not in _RETRY_STATUS:
                    raise UpstreamError(f"{host} {last_err}")
                if retry_after := resp.headers.get("retry-after"):
                    try:
                        await asyncio.sleep(min(float(retry_after), 30))
                        continue
                    except ValueError:
                        pass
            delay = min(2**attempt + random.uniform(0, 1), 20)
            log(f"{host} attempt {attempt + 1}/{attempts} failed ({last_err[:120]}); retrying in {delay:.1f}s")
            await asyncio.sleep(delay)
        raise UpstreamError(f"{host} unavailable after {attempts} attempts ({last_err[:200]})")

    async def get_json(self, url: str, **kw: Any) -> Any:
        return json.loads(await self.request("GET", url, **kw))

    async def post_json(self, url: str, body: Any, **kw: Any) -> Any:
        return json.loads(await self.request("POST", url, json_body=body, **kw))

    async def aclose(self) -> None:
        await self._client.aclose()


_http: Http | None = None


def http() -> Http:
    global _http
    if _http is None:
        _http = Http()
    return _http
