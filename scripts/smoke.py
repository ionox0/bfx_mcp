"""Call every tool through a real MCP client session (in-memory transport).

    uv run python scripts/smoke.py [variant ...]
"""

import asyncio
import json
import sys
import time

from mcp import Client

from genomics_mcp.server import mcp


def show(title: str, result) -> None:
    body = result.structured_content if getattr(result, "structured_content", None) else [c.text for c in result.content]
    text = json.dumps(body, indent=1, default=str)
    flag = " (isError)" if result.is_error else ""
    print(f"\n=== {title}{flag}\n{text[:1800]}{' ...' if len(text) > 1800 else ''}")


async def main(variants: list[str]) -> None:
    async with Client(mcp) as client:
        tools = await client.list_tools()
        print("tools:", [t.name for t in tools.tools])
        print("resources:", [r.uri for r in (await client.list_resources()).resources])
        print("templates:", [t.uri_template for t in (await client.list_resource_templates()).resource_templates])
        print("prompts:", [p.name for p in (await client.list_prompts()).prompts])

        show("get_gene_constraint BRCA1", await client.call_tool("get_gene_constraint", {"gene_symbol": "BRCA1"}))
        for v in variants:
            for tool in ("normalize_variant", "get_population_frequency", "get_clinvar_classification", "annotate_variant"):
                t0 = time.monotonic()
                res = await client.call_tool(tool, {"variant": v})
                show(f"{tool}({v}) {time.monotonic() - t0:.1f}s", res)


if __name__ == "__main__":
    asyncio.run(main(sys.argv[1:] or ["NM_000492.4:c.1521_1523del", "rs80357258", "17-43045712-T-C", "1-1-A-T"]))
