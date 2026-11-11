# -*- coding: utf-8 -*-
"""Open Targets GraphQL 查询小工具（避免 shell 引号问题）。"""
import json
import sys
import time
from urllib.request import Request, urlopen

URL = 'https://api.platform.opentargets.org/api/v4/graphql'


def query(q, timeout=90, retries=3):
    body = json.dumps({'query': q}).encode('utf-8')
    last = None
    for k in range(retries):
        try:
            req = Request(URL, data=body,
                          headers={'Content-Type': 'application/json',
                                   'User-Agent': 'A1-MR/1.0'})
            with urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception as e:
            last = e
            time.sleep(2 * (k + 1))
    raise last


if __name__ == '__main__':
    q = sys.stdin.read()
    print(json.dumps(query(q), ensure_ascii=False, indent=2))
