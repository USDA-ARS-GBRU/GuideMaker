# Edition

Product documentation template for Jekyll. Browse through a [live demo](https://long-pig.cloudvent.net/).
Start documenting your product, application, service or website with this configurable theme.

![Edition template screenshot](images/_screenshot.png)

<<<<<<< HEAD
Edition was made by [CloudCannon](http://cloudcannon.com/), the Cloud CMS for Jekyll.
=======
>>>>>>> 16b5fff1ffb3fb4415509a8b353841d1dc428997

Find more templates, themes and step-by-step Jekyll tutorials at [CloudCannon Academy](https://learn.cloudcannon.com/).

<<<<<<< HEAD
## Features
=======
**Methods to access GuideMaker**


1.[Web Application](https://guidemaker.app.scinet.usda.gov)      |2.[CyCverse Discovery Environment](https://cyverse.org/discovery-environment)
:---------------------:|:-------------------------:
[![Image of Guidemaker Web App](https://raw.githubusercontent.com/USDA-ARS-GBRU/GuideMaker/main/guidemaker/data/scinet.png)](https://guidemaker.app.scinet.usda.gov) | [![Image of Guidemaker Web App](https://raw.githubusercontent.com/USDA-ARS-GBRU/GuideMaker/main/guidemaker/data/cyverse.png)](https://cyverse.org/discovery-environment)


>>>>>>> 16b5fff1ffb3fb4415509a8b353841d1dc428997

* Two column layout
* Full text search
* Pre-styled components
* Auto-generated navigation based on category
* Optimised for editing in [CloudCannon](http://cloudcannon.com/)
* Change log
* RSS/Atom feed
* SEO tags
* Google Analytics

## Setup

<<<<<<< HEAD
1. Add your site and author details in `_config.yml`.
2. Get a workflow going to see your site's output (with [CloudCannon](https://app.cloudcannon.com/) or Jekyll locally).
=======
### 3.Command Line
>>>>>>> 16b5fff1ffb3fb4415509a8b353841d1dc428997

## Develop

Edition was built with [Jekyll](http://jekyllrb.com/) version 3.3.1, but should support newer versions as well.

Install the dependencies with [Bundler](http://bundler.io/):

~~~bash
$ bundle install
~~~

Run `jekyll` commands through Bundler to ensure you're using the right versions:

~~~bash
$ bundle exec jekyll serve
~~~

## Editing

Edition is already optimised for adding, updating and removing documentation pages in CloudCannon.

### Documentation pages

* Add, update or remove a documentation page in the *Documentation* collection.
* Change the category of a documentation page to move it to another section in the navigation.
* Documentation pages are organised in the navigation by category, with URLs based on the path inside the `_docs` folder.

### Change log

* Add, update or remove change log entries from your posts.
* Tag entries as minor or major in the front matter.

### Search

* Add `excluded_in_search: true` to any documentation page's front matter to exclude that page in the search results.

### Navigation

* Change `site.show_full_navigation` to control all or only the current navigation group being open.
