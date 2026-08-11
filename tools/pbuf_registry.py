#!/usr/bin/env python3
"""Query, validate, and render the permanent PBUF scientific memory."""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pbuf.registry.search import load, search
from pbuf.registry.validate import validate
from pbuf.registry.render import render_registry, render_timeline

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'docs/PBUF_MECHANISM_REGISTRY.json'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['search','target','attempt','history','canonical','equivalents','reopen','validate','render'])
    parser.add_argument('query', nargs='*'); parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args(); data = load(REGISTRY); query = ' '.join(args.query)
    if args.command == 'validate':
        errors = validate(data); print(json.dumps({'valid': not errors, 'errors': errors}, indent=2)); return 1 if errors else 0
    if args.command == 'render':
        commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
        date = subprocess.check_output(['git','show','-s','--format=%as','HEAD'], cwd=ROOT, text=True).strip()
        (ROOT/'docs/PBUF_MECHANISM_REGISTRY.md').write_text(render_registry(data, commit, date))
        (ROOT/'docs/PBUF_DEVELOPMENT_TIMELINE.md').write_text(render_timeline(data)); print('rendered'); return 0
    if args.command == 'target': ts = [t for t in data['targets'] if t['target_id'] == query]
    elif args.command == 'attempt':
        ts = []; ats = [a for a in data['attempts'] if query.lower() in hay(a)]
        print_records(ts, ats, args.verbose); return 0
    elif args.command == 'equivalents':
        normalized = query.lower().replace(' ', '_')
        matching_ids = {a['attempt_id'] for a in data['attempts'] if query.lower() in json.dumps(a).lower() or normalized in json.dumps(a).lower()}
        rel = [e for e in data.get('equivalences', []) if query.lower() in json.dumps(e).lower() or normalized in json.dumps(e).lower() or e['source'] in matching_ids or e['target'] in matching_ids]
        print(json.dumps(rel, indent=2)); return 0
    elif args.command == 'canonical': ts = [t for t in data['targets'] if t['current_status']=='CANONICAL' and query.lower() in json.dumps(t).lower()]
    elif args.command == 'reopen': ts = [t for t in data['targets'] if query.lower() in json.dumps(t).lower()]
    else: ts, ats = search(data, query); print_records(ts, ats, args.verbose); return 0
    ats = [a for a in data['attempts'] if a['target_id'] in {t['target_id'] for t in ts}]; print_records(ts, ats, args.verbose)

def hay(x): return json.dumps(x).lower()
def print_records(targets, attempts, verbose):
    by_target = {t['target_id']: t for t in targets}
    for t in targets:
        print(f"TARGET: {t['canonical_name']}\nSTATUS: {t['current_status']}\n")
    for a in attempts:
        print(f"{a.get('dev') or a.get('pr') or 'Pre-ledger'}  {a['name']}  {a['result']} / {a['current_status']}")
        if verbose: print(json.dumps(a, indent=2))
if __name__ == '__main__': sys.exit(main())
