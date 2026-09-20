import json
import urllib.request
import ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=45, context=ctx) as r:
        return r.read().decode("utf-8", "replace")

for url in [
    "https://www.ebi.ac.uk/gwas/rest/api/publications/33729212/studies",
    "https://www.ebi.ac.uk/gwas/rest/api/studies?q=33729212&size=10",
]:
    print(">>>", url)
    try:
        d = json.loads(get(url))
        recs = d.get("_embedded", {}).get("studies", [])
        if not recs and isinstance(d, list):
            recs = d
        for s in recs:
            anc = s.get("ancestries", [])
            eur = any("european" in ((a.get("type") or "") + str(a.get("ancestralGroups") or [])).lower() for a in anc)
            ss = any(a.get("summaryStats") for a in anc)
            print(f"  {s.get('accessionId'):12s} | trait={str(s.get('traitName',''))[:40]:40s} | "
                  f"EUR={eur} ss={ss} | n={str(s.get('initialSampleSize',''))[:60]}")
    except Exception as e:
        print("  ERR", repr(e))
    print()
