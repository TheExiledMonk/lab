#!/usr/bin/env python3
import argparse,functools,http.server,json,socketserver
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument("dataset",type=Path); p.add_argument("--port",type=int,default=0); a=p.parse_args()
root=Path(__file__).resolve().parent; dataset=a.dataset.resolve()
class Handler(http.server.SimpleHTTPRequestHandler):
    def translate_path(self,path):
        if path.startswith("/data/"): return str(dataset/path.removeprefix("/data/"))
        return str(root/(path.lstrip("/") or "index.html"))
with socketserver.TCPServer(("127.0.0.1",a.port),Handler) as s:
    print(f"http://127.0.0.1:{s.server_address[1]}/",flush=True); s.serve_forever()
