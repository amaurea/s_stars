all: table.txt table.html table.wiki

table.txt: params.txt reformat
	./reformat -f ascii $< $@

table.html: params.txt reformat
	./reformat -f html $< $@

table.wiki: params.txt reformat
	./reformat -f wiki $< $@

reformat: reformat.c
	gcc -g -o $@ $< -lm

clean:
	rm -f table.txt table.html table.wiki
