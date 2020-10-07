
SHELL := /bin/zsh

## Compilators
##============
TEX=pdflatex
BIB=bibtex
PYTHON3=python3

## File to build
##==============
SRC=main.tex
SRCBIB=bibliography.bib
TARGET=$(SRC:.tex=.pdf)

## Build log files
##================
SRCLOG=$(SRC:.tex=.log)
SRCAUX=$(SRC:.tex=.aux)
SRCIDX=$(SRC:.tex=.idx)
SRCILG=$(SRC:.tex=.ilg)
SRCIND=$(SRC:.tex=.ind)
SRCTOC=$(SRC:.tex=.toc)
SRCBBL=$(SRC:.tex=.bbl)
SRCBLG=$(SRC:.tex=.blg)
SRCNAV=$(SRC:.tex=.nav)
SRCOUT=$(SRC:.tex=.out)
SRCSNM=$(SRC:.tex=.snm)
SRCBLX=$(SRC:.tex=-blx.bib)
SRCXML=$(SRC:.tex=.run.xml)


## Figures
##========
FIGPY=$(wildcard figures/*.py)
FIGPNG=$(FIGPY:%.py=%.png)


## Make rules
##===========

all: $(TARGET)

$(TARGET): $(SRC) $(SRCBIB) $(FIGPNG)
	@$(TEX) $(SRC)
	@$(BIB) $(SRCAUX)
	@$(TEX) $(SRC)
	@$(TEX) $(SRC)
	@cp $(TARGET) ../doc.pdf
	@mv $(SRCLOG) logs 2> /dev/null | true
	@mv $(SRCAUX) logs 2> /dev/null | true
	@mv $(SRCIDX) logs 2> /dev/null | true
	@mv $(SRCILG) logs 2> /dev/null | true
	@mv $(SRCIND) logs 2> /dev/null | true
	@mv $(SRCTOC) logs 2> /dev/null | true
	@mv $(SRCBBL) logs 2> /dev/null | true
	@mv $(SRCBLG) logs 2> /dev/null | true
	@mv $(SRCNAV) logs 2> /dev/null | true
	@mv $(SRCOUT) logs 2> /dev/null | true
	@mv $(SRCSNM) logs 2> /dev/null | true
	@mv $(SRCBLX) logs 2> /dev/null | true
	@mv $(SRCXML) logs 2> /dev/null | true


## Figures rules
##==============

figures/%.png: figures/%.py
	@$(PYTHON3) $<
#	@echo compiling "$< ==> $@"


## Clean rules
##============

clean:
	@rm -f $(TARGET) 2> /dev/null | true

clean-figures:
	@rm -f $(FIGPNG) 2> /dev/null | true

mrproper:
	@rm -f $(SRCLOG) 2> /dev/null | true
	@rm -f $(SRCAUX) 2> /dev/null | true
	@rm -f $(SRCIDX) 2> /dev/null | true
	@rm -f $(SRCILG) 2> /dev/null | true
	@rm -f $(SRCIND) 2> /dev/null | true
	@rm -f $(SRCTOC) 2> /dev/null | true
	@rm -f $(SRCBBL) 2> /dev/null | true
	@rm -f $(SRCBLG) 2> /dev/null | true
	@rm -f $(SRCNAV) 2> /dev/null | true
	@rm -f $(SRCOUT) 2> /dev/null | true
	@rm -f $(SRCSNM) 2> /dev/null | true
	@rm -f $(SRCBLX) 2> /dev/null | true
	@rm -f $(SRCXML) 2> /dev/null | true
	@rm -f logs/$(SRCLOG) 2> /dev/null | true
	@rm -f logs/$(SRCAUX) 2> /dev/null | true
	@rm -f logs/$(SRCIDX) 2> /dev/null | true
	@rm -f logs/$(SRCILG) 2> /dev/null | true
	@rm -f logs/$(SRCIND) 2> /dev/null | true
	@rm -f logs/$(SRCTOC) 2> /dev/null | true
	@rm -f logs/$(SRCBBL) 2> /dev/null | true
	@rm -f logs/$(SRCBLG) 2> /dev/null | true
	@rm -f logs/$(SRCNAV) 2> /dev/null | true
	@rm -f logs/$(SRCOUT) 2> /dev/null | true
	@rm -f logs/$(SRCSNM) 2> /dev/null | true
	@rm -f logs/$(SRCBLX) 2> /dev/null | true
	@rm -f logs/$(SRCXML) 2> /dev/null | true
	@rm -f $(FIGPNG) 2> /dev/null | true
	@rm -f $(TARGET) 2> /dev/null | true

