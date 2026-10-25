"""
A teaching MCP server: a tiny "notebook" that shows off every core MCP idea.

READ THIS FIRST
===============
MCP (Model Context Protocol) is a standard way for an AI app (the *host*, e.g.
Claude Code or Claude Desktop) to plug in outside capabilities. The host starts
an MCP *client* for each server it talks to, and the two exchange JSON-RPC 2.0
messages.

A server can offer three kinds of things ("primitives"):

  1. TOOLS      -> functions the MODEL decides to call      (like POST: do something)
  2. RESOURCES  -> data the APP/user attaches as context    (like GET: read something)
  3. PROMPTS    -> reusable templates the USER picks        (like slash commands)

This file defines one or more of each. The SDK (`MCPServer`) turns your plain
Python functions into protocol messages:
  * the function name          -> the tool/prompt name
  * the docstring              -> the description the model reads
  * the type hints             -> a JSON Schema for the arguments
  * the return value           -> the result content

IMPORTANT stdio RULE: with the stdio transport, stdout IS the protocol channel.
Never print() to stdout in a server - it corrupts the JSON-RPC stream. Log to
stderr instead (that's what `log()` below does).

Run `uv run python raw_client.py` to watch the actual messages go back and forth.
"""

import sys
from datetime import datetime

from mcp.server.mcpserver import Context, MCPServer
from mcp.server.mcpserver.exceptions import ToolError


def log(msg: str) -> None:
    """Debug output goes to STDERR, never stdout (see the stdio rule above)."""
    print(f"[server] {msg}", file=sys.stderr, flush=True)


# The server object. `name` and `instructions` are sent to the client during the
# `initialize` handshake; hosts often show `instructions` to the model.
mcp = MCPServer(
    name="learn-mcp",
    instructions=(
        "A small notebook. Use add_note to save notes, search_notes to find "
        "them, and read notes://all for everything."
    ),
)

# Our "database". A real server would talk to a DB, an API, the filesystem...
NOTES: dict[int, dict] = {
    1: {"title": "What is MCP?", "body": "A protocol connecting AI apps to tools and data.", "created": "2026-09-30T09:00:00"},
    2: {"title": "Transports", "body": "stdio for local processes, Streamable HTTP for remote servers.", "created": "2026-09-30T09:05:00"},
}


# ---------------------------------------------------------------------------
# 1. TOOLS - the model sees the name + docstring + argument schema, and decides
#    on its own when to call them. They can have side effects.
#    Wire messages: `tools/list` (discover) and `tools/call` (invoke).
# ---------------------------------------------------------------------------

@mcp.tool()
def add_note(title: str, body: str) -> str:
    """Save a new note to the notebook. Returns the new note's id."""
    note_id = max(NOTES, default=0) + 1
    NOTES[note_id] = {"title": title, "body": body, "created": datetime.now().isoformat(timespec="seconds")}
    log(f"add_note -> created note {note_id}")
    return f"Saved note {note_id}: {title!r}"


@mcp.tool()
def search_notes(query: str, limit: int = 5) -> list[dict]:
    """Find notes whose title or body contains `query` (case-insensitive)."""
    # Returning structured data (a list of dicts) makes the SDK send both a text
    # version and `structuredContent` that clients can parse reliably.
    q = query.lower()
    hits = [
        {"id": i, "title": n["title"], "body": n["body"]}
        for i, n in NOTES.items()
        if q in n["title"].lower() or q in n["body"].lower()
    ]
    log(f"search_notes({query!r}) -> {len(hits)} hit(s)")
    return hits[:limit]


@mcp.tool()
def delete_note(note_id: int) -> str:
    """Delete a note by id."""
    # Raising ToolError is how you report an *expected* failure. The SDK turns
    # it into a normal result with `isError: true` and your message, so the
    # MODEL sees it and can react (e.g. retry with a valid id). Any other
    # exception is treated as a crash and the model only sees a generic
    # "Error executing tool" - so your helpful message would be lost.
    if note_id not in NOTES:
        raise ToolError(f"No note with id {note_id}. Existing ids: {sorted(NOTES)}")
    title = NOTES.pop(note_id)["title"]
    return f"Deleted note {note_id}: {title!r}"


@mcp.tool()
async def count_slowly(to: int, ctx: Context) -> str:
    """Count up to `to`, reporting progress along the way (demo of Context)."""
    # Asking for a `Context` parameter gives you a handle back to the client.
    # It is hidden from the model's argument schema. Here we use it to send
    # `notifications/progress` messages while the tool is still running - but
    # only if the client asked for progress by sending a `progressToken`.
    import asyncio

    for i in range(1, to + 1):
        await ctx.report_progress(progress=i, total=to, message=f"at {i}")
        await asyncio.sleep(0.2)
    return f"Counted to {to}"


# ---------------------------------------------------------------------------
# 2. RESOURCES - read-only data identified by a URI. The *application* (or the
#    user, e.g. typing @ in Claude Code) chooses which to load into context.
#    Wire messages: `resources/list`, `resources/templates/list`, `resources/read`.
# ---------------------------------------------------------------------------

@mcp.resource("notes://all", mime_type="text/markdown")
def all_notes() -> str:
    """Every note in the notebook, as markdown."""
    return "\n\n".join(f"## [{i}] {n['title']}\n{n['body']}" for i, n in NOTES.items())


# A URI with {placeholders} becomes a *resource template*: one definition that
# covers many URIs (notes://1, notes://2, ...). The placeholder maps to the arg.
@mcp.resource("notes://{note_id}")
def one_note(note_id: str) -> str:
    """A single note by id."""
    n = NOTES.get(int(note_id))
    return f"# {n['title']}\n\n{n['body']}\n\n_created {n['created']}_" if n else "Not found"


# ---------------------------------------------------------------------------
# 3. PROMPTS - parameterised message templates that the USER explicitly picks
#    (Claude Code shows them as /mcp__learn-mcp__<name> slash commands).
#    Wire messages: `prompts/list`, `prompts/get`.
# ---------------------------------------------------------------------------

@mcp.prompt()
def summarize_notes(style: str = "bullet points") -> str:
    """Ask the model to summarize the whole notebook."""
    # A prompt only *produces messages*; it doesn't call the model itself.
    # The host puts these messages into the conversation.
    return f"Summarize these notes as {style}:\n\n{all_notes()}"


if __name__ == "__main__":
    log("starting on stdio - waiting for a client to send `initialize`")
    # stdio transport: the client launches this process and talks over
    # stdin/stdout. For a remote server you'd use mcp.run("streamable-http").
    mcp.run("stdio")
