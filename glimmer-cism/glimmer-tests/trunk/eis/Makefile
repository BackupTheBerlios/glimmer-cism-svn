SUBDIRS = EISMINT-1 EISMINT-2 eis

subdirs:	$(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@

clean:
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir clean; \
		done

.PHONY:		clean subdirs $(SUBDIRS)
