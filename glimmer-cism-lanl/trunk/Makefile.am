SUBDIRS = m4macros src 
EXTRA_DIST = bootstrap tests2
bin_SCRIPTS=glimmer-config

# Make sure we never package up '.svn' directories
dist-hook:
	find $(distdir) -name '.svn' | xargs rm -rf

