EXTRA_DIST = test_transport.sh trans_velo.config.in trans_velo.py

BUILT_SOURCES = trans_velo.eps

trans_velo.ps:	test_transport.sh trans_velo.config.in trans_velo.py ../src/test_transport
		bash ./test_transport.sh

include $(top_srcdir)/doc/make_doc.rules
