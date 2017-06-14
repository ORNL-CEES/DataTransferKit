# Inspired from https://gist.github.com/kakawait/9215487 visited on June 14th, 2017
# See https://github.com/sphinx-doc/sphinx/issues/1717#issuecomment-74420054
# and https://github.com/sphinx-doc/sphinx/issues/1717#issuecomment-74679763
#
# Sphinx extension to conditionaly inculde files

import os, re
from sphinx import addnodes

docs_to_remove = []

def setup(app):
    app.ignore = []
    app.connect('builder-inited', builder_inited)
    app.connect('env-get-outdated', env_get_outdated)
    app.connect('doctree-read', doctree_read)
    app.add_config_value("scope_donotinclude", [], '')

def builder_inited(app):
    for doc in app.env.found_docs:
        first_directive = None
        for suffix in app.env.config.source_suffix:
            _file = app.env.srcdir + os.sep + doc + suffix
            if os.path.isfile(_file):
                with open(_file, 'r') as fin:
                    first_directive = fin.readline() + fin.readline()
                if first_directive:
                    # check whether the first two lines of a file match
                    # ```
                    # .. meta::
                    #     :scope: <TAG>
                    # ```
                    # and add file if <TAG> belongs to the "donotinclude" list
                    m = re.match(r'^\.\. meta::\s+:scope: ([a-zA-Z0-9_-]+)', first_directive)
                    if m and m.group(1) in app.config.scope_donotinclude:
                        docs_to_remove.append(doc)
    app.env.found_docs.difference_update(docs_to_remove)

def env_get_outdated(app, env, added, changed, removed):
    added.difference_update(docs_to_remove)
    changed.difference_update(docs_to_remove)
    removed.update(docs_to_remove)
    return []

def doctree_read(app, doctree):
    for toctreenode in doctree.traverse(addnodes.toctree):
        for e in toctreenode['entries']:
            ref = str(e[1])
            if ref in docs_to_remove:
                toctreenode['entries'].remove(e)
