#
#  Makefile: generating HTML and PDF from RST
#
#  NB: not a good solution, just an initial hack.


SRC = ATSNativeSpec_dev.rst.in

RST = $(SRC:%.rst.in=%.rst)

# generic rules for targets

all : rst

%.rst:
	@./generate_ATSNativeSpec.py $< $@

rst : $(RST)

clean :
	-/bin/rm -f $(RST)
	-/bin/rm -f *~

