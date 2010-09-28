SUBDIRS = m4macros src tests/util
EXTRA_DIST = bootstrap 
bin_SCRIPTS=glimmer-config

# Add a target to run an automatic test
test::
	${prefix}/tests/util/automatic_test.sh ${prefix} ${FC} ${expt_type} ${expt_num}

# Make sure we never package up '.svn' directories
dist-hook:
	find $(distdir) -name '.svn' | xargs rm -rf

