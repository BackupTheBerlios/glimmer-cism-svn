SUBDIRS = EISMINT-1 EISMINT-2 lin-non_lin eis

subdirs:	$(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@
plots:
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir plots; \
		done

clean:
		for dir in $(SUBDIRS); do \
		  $(MAKE) -C $$dir clean; \
		done

.PHONY:		clean subdirs $(SUBDIRS)
