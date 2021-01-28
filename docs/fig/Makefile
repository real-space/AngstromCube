# Grace plots (xmgrace *.agr)
AGRFILES = $(wildcard *.agr)

# Fig drawings (xfig *.fig)
FIGFILES = $(wildcard *.fig)

TARGET = $(FIGFILES:%.fig=%.pdf) $(FIGFILES:%.fig=%.svg) $(AGRFILES:%.agr=%.pdf)

all: $(TARGET)

veryall: all

.PHONY: clean veryclean

%.eps:	%.fig
	fig2dev -L eps $< $@

%.pdf:	%.fig
	fig2dev -L pdf $< $@

%.svg:	%.fig
	fig2dev -L svg $< $@

%.png:	%.fig
	fig2dev -L png $< $@

%.pdf:	%.dot
	dot -Tpdf -o $@ $<

%.pdf:	%.gif
	convert $*.gif $*.pdf

%.eps:  %.agr
	gracebat -hdevice EPS $*.agr

%.pdf:  %.agr
	gracebat -hdevice EPS $*.agr && \
	epstopdf $*.eps && \
	rm -f $*.eps

clean:
	-rm -f *.bak

veryclean: clean
	-rm -f $(TARGET)
