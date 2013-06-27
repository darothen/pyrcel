import urllib2, urllib
import sys
from ghost import Ghost

import lxml.html
from lxml.cssselect import CSSSelector

MASTER_NODE = "legion.mit.edu"
API_URL = "http://%s:8989/job.html" % MASTER_NODE

def call_api(job_name):
    url = API_URL + "?name=%s" % job_name
    req = urllib2.Request(url)
    return urllib2.urlopen(req)

if __name__ == "__main__":

    job = "RunParcelModels@554:7c051:3b4b6"

    ghost = Ghost()
    url = API_URL + "?name=%s" % job
    ghost.open(url)
    ghost.wait_for_page_loaded()

    page_html = str(ghost.content)
    page_lines = page_html.split("\n")
    for line in page_lines:
        line = line.strip()
        if line.startswith("<div class=\"events\">"): break

    events = reversed(line.split("<div class=\"event\">")[1:])

    for event in events:
        divs = lxml.html.fromstring(event)
        ts_sel = CSSSelector("div.tstamp")
        node_sel = CSSSelector("div.node")
        text_sel = CSSSelector("div.text ")

        timestamp = ts_sel(divs)[0].text
        node =  node_sel(divs)[0].text

        text_elem = text_sel(divs)[0]
        text = text_elem.cssselect('pre')[0].text

        print "@ %s, %s >> %s" % (timestamp, node, text)





