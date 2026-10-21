import json, re, pathlib
d = pathlib.Path(__file__).parent
t = (d / "template.html").read_text()
data = (d / "data.json").read_text()
page = t.replace("/*__DATA__*/", data)
(d / "baltimore_heat_atlas.html").write_text(page)            # artifact body (skeleton added at publish)
head, _, body = page.partition("</style>")
site = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
        + head + "</style>\n</head>\n<body>\n" + body + "\n</body>\n</html>\n")
(d.parent / "docs").mkdir(exist_ok=True)
(d.parent / "docs" / "index.html").write_text(site)
print(len(page)//1024, "KB")
