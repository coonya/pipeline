import urllib2

address = 'http://admin-apps.webofknowledge.com/JCR/JCR?RQ=LIST_SUMMARY_JOURNAL&cursor=1781'
website = urllib2.urlopen(address)
website_html = website.read()
print website_html





print 'aaaaaaaaaaaaa'
