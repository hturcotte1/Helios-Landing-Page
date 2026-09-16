"""Render results.md (+ README) into a single self-contained HTML manuscript with MathJax (cdnjs)."""
import markdown, html, os, sys, re, datetime
root = os.path.join(os.path.dirname(__file__), "..")
md = open(os.path.join(root, "results.md")).read()
readme = open(os.path.join(root, "README.md")).read()
# protect math from the markdown converter: replace $$...$$ and $...$ with placeholders
holders = []
def keep(m):
    holders.append(m.group(0)); return f"\x00{len(holders)-1}\x00"
md2 = re.sub(r"\$\$.*?\$\$", keep, md, flags=re.S)
md2 = re.sub(r"(?<!\\)\$(?!\s)(.+?)(?<!\s)\$", keep, md2)
body = markdown.markdown(md2, extensions=["tables", "fenced_code"])
body = re.sub(r"\x00(\d+)\x00", lambda m: html.escape(holders[int(m.group(1))], quote=False).replace("&amp;", "&"), body)
readme_html = markdown.markdown(readme, extensions=["tables", "fenced_code"])
date = datetime.date.today().isoformat()
page = f"""<title>SYT Tree Automorphisms and Census</title>
<style>
:root{{--bg:#fbfaf7;--fg:#1c1b19;--muted:#5d5a53;--rule:#d9d5cc;--code:#f0ede6;--accent:#3b5b8a}}
@media (prefers-color-scheme: dark){{:root:not([data-theme="light"]){{--bg:#171614;--fg:#e8e4dc;--muted:#a9a49a;--rule:#3a3733;--code:#25231f;--accent:#8fb3e8}}}}
:root[data-theme="dark"]{{--bg:#171614;--fg:#e8e4dc;--muted:#a9a49a;--rule:#3a3733;--code:#25231f;--accent:#8fb3e8}}
body{{background:var(--bg);color:var(--fg);font:16px/1.55 Georgia,'Times New Roman',serif;margin:0;padding-block:2rem;padding-inline:max(16px,calc(50vw - 24rem))}}
h1,h2,h3{{line-height:1.2;font-family:Georgia,serif}} h1{{font-size:1.9rem}} h2{{font-size:1.4rem;margin-top:2.2rem;border-top:1px solid var(--rule);padding-top:1rem}} h3{{font-size:1.15rem;margin-top:1.6rem}}
p{{margin:.7rem 0}} a{{color:var(--accent)}} code{{background:var(--code);padding:.05em .3em;border-radius:3px;font-size:.9em}}
pre{{background:var(--code);padding:.8rem;overflow-x:auto;border-radius:4px}} pre code{{background:none;padding:0}}
table{{border-collapse:collapse;display:block;overflow-x:auto;max-width:100%}} td,th{{border:1px solid var(--rule);padding:.25rem .5rem}}
.meta{{color:var(--muted);font-size:.95rem}} .readme{{border:1px solid var(--rule);padding:1rem 1.2rem;border-radius:6px;margin-bottom:2rem;background:var(--code)}}
mjx-container{{overflow-x:auto;overflow-y:hidden;max-width:100%}}
</style>
<script>window.MathJax={{tex:{{inlineMath:[['$','$']],displayMath:[['$$','$$']]}},options:{{skipHtmlTags:['script','noscript','style','textarea','pre','code']}}}};</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-chtml.js"></script>
<p class="meta">Manuscript rendered {date} from <code>syt_tree/results.md</code> (branch <code>claude/syt-tree-automorphisms-zvgoes</code>). The README summary precedes the full write-up.</p>
<div class="readme">{readme_html}</div>
{body}
"""
out = os.path.join(root, "manuscript.html")
open(out, "w").write(page)
print("wrote", out, len(page), "bytes")
