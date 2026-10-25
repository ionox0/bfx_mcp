# learn-mcp

A small, heavily commented MCP server (a "notebook") for learning how the
Model Context Protocol works, plus a raw client that shows every message.

| File | What it teaches |
|---|---|
| `server.py` | How to **build** a server: tools, resources, resource templates, prompts, errors, progress |
| `raw_client.py` | What actually goes **over the wire**. Uses no MCP library, only `subprocess` + `json` |

## Quick start

```bash
uv run python raw_client.py           # watch a full session scroll by
uv run python raw_client.py --step    # pause after each step
```

## The mental model

```
┌──────────── Host (Claude Code, Claude Desktop, …) ─────────────┐
│   LLM  ◄──►  MCP client ─────────┐                             │
└──────────────────────────────────┼─────────────────────────────┘
                                   │  JSON-RPC 2.0 messages
                                   │  (stdio pipe, or HTTP)
                            ┌──────▼──────┐
                            │ MCP server  │  ← server.py
                            │ tools       │
                            │ resources   │
                            │ prompts     │
                            └─────────────┘
```

* The **host** is the AI app. It runs one **client** per server.
* The **server** exposes capabilities. It never talks to the LLM directly; the
  host decides what the model sees.
* **Transport**: `stdio` means the host launches the server as a subprocess and
  exchanges one JSON message per line on stdin and stdout. For remote servers,
  **Streamable HTTP** carries the same messages over HTTP.

### The three primitives

| Primitive | Who decides to use it | Analogy | Methods |
|---|---|---|---|
| **Tools** | the **model** | POST, an action | `tools/list`, `tools/call` |
| **Resources** | the **app or user** | GET, read-only data by URI | `resources/list`, `resources/templates/list`, `resources/read` |
| **Prompts** | the **user** | slash-command templates | `prompts/list`, `prompts/get` |

### The lifecycle

1. `initialize` (request): client and server agree on a protocol version and
   exchange capabilities.
2. `notifications/initialized` (notification): the handshake is done.
3. Normal traffic: `*/list` to discover capabilities, then `tools/call`,
   `resources/read` and the rest.
4. Shutdown: for stdio, the client closes stdin.

### Message types (JSON-RPC 2.0)

* **Request**: has an `id` and expects exactly one response with the same `id`.
* **Response**: contains either `result` or `error`.
* **Notification**: has no `id`, gets no reply, and can go in either direction
  (for example `notifications/progress`).

### Two kinds of errors (step 4 vs. step 8 in the client)

* **Tool error**: a normal `result` with `isError: true`. The model reads the
  message and can recover. In SDK 2.x, raise `ToolError` for this. Any other
  exception is hidden behind a generic "Error executing tool".
* **Protocol error**: a JSON-RPC `error` object with a code (for example
  `-32601 Method not found`). It means the request itself was invalid.

### Gotchas

* **Never `print()` to stdout in a stdio server.** stdout carries the protocol,
  so log to stderr.
* The **docstring and type hints are the tool's interface to the model**. Write
  them for an LLM reader.

## Use it from Claude Code

```bash
claude mcp add learn-mcp -- uv run --directory /Users/ianjohnson/Desktop/mcp_server python server.py
```

Then, inside `claude`:
* `/mcp`: check that it is connected and see its tools
* "save a note that MCP is neat": the model calls `add_note`
* type `@` and pick `learn-mcp:notes://all`: this attaches a **resource**
* `/mcp__learn-mcp__summarize_notes`: this runs the **prompt**

To remove it: `claude mcp remove learn-mcp`.

## Other ways to poke at it

* **MCP Inspector** (a web UI, needs Node): `npx @modelcontextprotocol/inspector uv run python server.py`
* Edit `server.py`: add a tool, rerun `raw_client.py`, and see how its schema appears in `tools/list`.

Built with the official Python SDK, `mcp` 2.2.0 (`MCPServer`, formerly `FastMCP`).
