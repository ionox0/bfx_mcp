"""
A bare-bones MCP client that uses NO MCP library - just subprocess + json.

It launches server.py, then walks through a full MCP session step by step,
printing every JSON-RPC message in both directions so you can see exactly what
a host like Claude Code does behind the scenes.

    uv run python raw_client.py          # run straight through
    uv run python raw_client.py --step   # pause after each step (press Enter)

Colour key:  --> client to server    <-- server to client
"""

import json
import subprocess
import sys
import threading

STEP = "--step" in sys.argv
BLUE, GREEN, YELLOW, DIM, BOLD, RESET = "\033[34m", "\033[32m", "\033[33m", "\033[2m", "\033[1m", "\033[0m"

# 1) The client LAUNCHES the server as a child process. This is the "stdio"
#    transport: we write JSON to its stdin, it writes JSON to its stdout, one
#    message per line. (Its stderr is just logs; we let it pass through.)
server = subprocess.Popen(
    [sys.executable, "server.py"],
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    text=True,
    bufsize=1,
)

next_id = 0


def show(arrow: str, colour: str, msg: dict) -> None:
    print(f"{colour}{arrow} {json.dumps(msg, indent=2)}{RESET}")


def send(msg: dict) -> None:
    show("-->", BLUE, msg)
    server.stdin.write(json.dumps(msg) + "\n")
    server.stdin.flush()


def request(method: str, params: dict | None = None) -> dict:
    """A REQUEST has an `id` and expects exactly one response with that id."""
    global next_id
    next_id += 1
    msg = {"jsonrpc": "2.0", "id": next_id, "method": method}
    if params is not None:
        msg["params"] = params
    send(msg)
    # Read lines until we get the response matching our id. Anything else that
    # arrives first is a NOTIFICATION from the server (no id, e.g. progress).
    while True:
        incoming = json.loads(server.stdout.readline())
        if incoming.get("id") == next_id:
            show("<--", GREEN, incoming)
            return incoming
        show("<-- (notification)", YELLOW, incoming)


def notify(method: str, params: dict | None = None) -> None:
    """A NOTIFICATION has no `id`: fire-and-forget, no response comes back."""
    msg = {"jsonrpc": "2.0", "method": method}
    if params is not None:
        msg["params"] = params
    send(msg)


def step(title: str, explanation: str) -> None:
    print(f"\n{BOLD}{'=' * 78}\n{title}{RESET}\n{DIM}{explanation}{RESET}\n")


def pause() -> None:
    if STEP:
        input(f"{DIM}[Enter to continue]{RESET}")


# ---------------------------------------------------------------------------
step(
    "STEP 1 - Handshake: initialize",
    "The first thing any client does. Both sides say which protocol version they\n"
    "speak and which optional features ('capabilities') they support. The server\n"
    "replies with its name, its capabilities (tools? resources? prompts?) and\n"
    "optional instructions for the model.",
)
init = request("initialize", {
    "protocolVersion": "2025-11-25",
    "capabilities": {},
    "clientInfo": {"name": "raw-client", "version": "0.1"},
})
print(f"{DIM}Server agreed on protocol {init['result']['protocolVersion']} and offers: "
      f"{', '.join(k for k in init['result']['capabilities'] if k != 'experimental')}{RESET}")
notify("notifications/initialized")
print(f"{DIM}^ 'initialized' notification: handshake done, normal traffic may begin.{RESET}")
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 2 - Discover tools: tools/list",
    "The host asks what tools exist. Each has a name, a description (from the\n"
    "Python docstring) and an inputSchema (JSON Schema, from the type hints).\n"
    "The host hands exactly this to the LLM - it's all the model knows about a tool.",
)
tools = request("tools/list")
print(f"{DIM}Tools the model would see: {[t['name'] for t in tools['result']['tools']]}{RESET}")
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 3 - Call a tool: tools/call",
    "When the MODEL decides to use a tool, the host sends tools/call with the\n"
    "arguments the model generated. The result's `content` goes back to the model.",
)
request("tools/call", {"name": "add_note", "arguments": {"title": "Hello", "body": "Written over raw JSON-RPC!"}})
request("tools/call", {"name": "search_notes", "arguments": {"query": "json"}})
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 4 - A tool error",
    "Note the difference: this is a *successful* JSON-RPC response whose result\n"
    "has isError: true. The error text is shown to the model so it can recover.\n"
    "(Protocol-level errors, like an unknown method, use a JSON-RPC `error` instead.)",
)
request("tools/call", {"name": "delete_note", "arguments": {"note_id": 99}})
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 5 - Progress notifications",
    "We attach a progressToken in _meta. While the tool runs, the server streams\n"
    "notifications/progress messages (yellow, no id) BEFORE the final response.",
)
request("tools/call", {"name": "count_slowly", "arguments": {"to": 3}, "_meta": {"progressToken": "count-1"}})
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 6 - Resources: resources/list, resources/templates/list, resources/read",
    "Resources are data addressed by URI. Fixed ones show up in resources/list;\n"
    "URI patterns like notes://{note_id} show up as templates. The APP (or you,\n"
    "via @-mentions in Claude Code) decides to read them - not the model.",
)
request("resources/list")
request("resources/templates/list")
request("resources/read", {"uri": "notes://3"})
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 7 - Prompts: prompts/list, prompts/get",
    "Prompts are templates the USER picks (slash commands in Claude Code).\n"
    "prompts/get returns ready-made messages; the host drops them into the chat.",
)
request("prompts/list")
request("prompts/get", {"name": "summarize_notes", "arguments": {"style": "a haiku"}})
pause()

# ---------------------------------------------------------------------------
step(
    "STEP 8 - A protocol-level error",
    "Calling a method the server doesn't know produces a JSON-RPC `error` object\n"
    "(with a standard code like -32601 'Method not found') instead of a `result`.",
)
request("does/not/exist")

# ---------------------------------------------------------------------------
step("DONE - Shutdown", "For stdio, shutting down is just closing stdin; the server exits.")
server.stdin.close()
server.wait(timeout=5)
print(f"Server exited with code {server.returncode}")
