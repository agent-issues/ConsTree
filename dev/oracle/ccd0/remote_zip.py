#!/usr/bin/env python3
"""List / extract members of a remote zip via HTTP range requests (stdlib only).
usage: remote_zip.py URL list [pattern]
       remote_zip.py URL get OUTDIR member [member ...]"""
import io, sys, os, fnmatch, zipfile, urllib.request, shutil

class HTTPRangeFile(io.RawIOBase):
    def __init__(self, url):
        self.orig = url
        self.resolve()
        self.pos = 0
    def resolve(self):
        # Signed redirect targets expire quickly; re-resolve on demand.
        req = urllib.request.Request(self.orig, headers={"Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            self.url = r.geturl()
            self.size = int(r.headers["Content-Range"].split("/")[-1])
            r.read()
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        if whence == 0: self.pos = off
        elif whence == 1: self.pos += off
        else: self.pos = self.size + off
        return self.pos
    def read(self, n=-1):
        if n is None or n < 0: n = self.size - self.pos
        if n == 0 or self.pos >= self.size: return b""
        end = min(self.size, self.pos + n) - 1
        for attempt in range(8):
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                break
            except Exception:
                if attempt == 7: raise
                try:
                    self.resolve()
                except Exception:
                    pass
        self.pos += len(data)
        return data
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)

url, mode = sys.argv[1], sys.argv[2]
f = io.BufferedReader(HTTPRangeFile(url), buffer_size=8 << 20)
z = zipfile.ZipFile(f)
if mode == "list":
    pat = sys.argv[3] if len(sys.argv) > 3 else "*"
    for i in z.infolist():
        if fnmatch.fnmatch(i.filename, pat):
            print(i.file_size, i.compress_size, i.filename)
else:
    out = sys.argv[3]
    for m in sys.argv[4:]:
        dst = os.path.join(out, m.replace("/", "__"))
        with z.open(m) as src, open(dst, "wb") as d:
            shutil.copyfileobj(src, d, 8 << 20)
        print("extracted", m, os.path.getsize(dst))
