#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals

AUTHOR = u'Andrew Mellor'
SITENAME = u'Empirical State of Mind'
HIDE_SITENAME = False
SITEURL = ''

PATH = 'content'

# Timezones
TIMEZONE = 'Europe/Paris'
DEFAULT_DATE_FORMAT = '%d %b, %Y'
DEFAULT_LANG = u'en'

# Set the article URL
ARTICLE_URL = 'blog/{date:%Y}/{date:%m}/{date:%d}/{slug}/'
ARTICLE_SAVE_AS = 'blog/{date:%Y}/{date:%m}/{date:%d}/{slug}/index.html'
#ARTICLE_URL = 'blog/{slug}.html'
#ARTICLE_SAVE_AS = 'blog/{slug}.html'

# Title menu options
MENUITEMS = []
NEWEST_FIRST_ARCHIVES = True
DISPLAY_CATEGORIES_ON_MENU = False
PAGE_ORDER_BY = 'order'

STATIC_PATHS = ['images', 'figures', 'code', 'notebooks', 'pdfs', 'favicon.png','extra/CNAME', 'data']
EXTRA_PATH_METADATA = {'extra/CNAME': {'path': 'CNAME'},}
CODE_DIR = 'code'
NOTEBOOK_DIR = 'notebooks'

# Theming + Plugins
THEME = "pelican-bootstrap3"
BOOTSTRAP_THEME = "flatly"
PLUGIN_PATHS = ['pelican-plugins']
PLUGINS = ['summary', 'liquid_tags.img', 'liquid_tags.video',
'liquid_tags.include_code', 'liquid_tags.notebook',
'liquid_tags.literal', 'render_math']

# Theme Extras
#SITELOGO = 'images/network.png'
#SITELOGO_SIZE 
FAVICON = 'images/network.png'

DEFAULT_PAGINATION = False

# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True

# Sidebar
SOCIAL = (('twitter', 'https://twitter.com/EmpiricalSOM'),
          ('linkedin', 'https://www.linkedin.com/pub/andrew-mellor/49/74b/842'),
          ('github', 'https://github.com/Andy1080'),)
DISPLAY_TAGS_ON_SIDEBAR = True
DISPLAY_CATEGORIES_ON_SIDEBAR = False
TAG_CLOUD_MAX_ITEMS = 10
DISPLAY_RECENT_POSTS_ON_SIDEBAR = True
RECENT_POST_COUNT = 5
HIDE_SIDEBAR = False


# Social - added with publishing
#DISQUS_SITENAME = 'empiricalstateofmind'
#ADDTHIS_PROFILE = 'ra-5477a96408ae534a'
#GOOGLE_ANALYTICS = 'UA-49010624-3'
#GOOGLE_ANALYTICS_UNIVERSAL = 'UA-49010624-3'
#GOOGLE_ANALYTICS_UNIVERSAL_PROPERTY = 'Empirical State of Mind'

#DISQUS_NO_ID = True
DISQUS_DISPLAY_COUNTS = True
USE_OPEN_GRAPH = True
TWITTER_CARDS = True
