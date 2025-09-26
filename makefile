build:
	rm -f ./a.out
	g++ -std=c++20 -g -Wall -pedantic -Werror main.cpp cleaner.o visualize.o -Wno-unused-variable -Wno-unused-function

run:
	./a.out

clean:
	rm -f ./a.out
	rm -f ecg_plot_*.png

valgrind:
	rm -f ./a.out
	g++ -std=c++20 -g -Wall -pedantic -Werror main.cpp cleaner.o visualize.o -Wno-unused-variable -Wno-unused-function
	valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./a.out

submit:
	/gradescope/gs submit course_id assignment_id files
  
commit:
	git add .
	git commit -m "$(msg)"

push:
	git push origin main

#
# Generate signal with optional rhythm and beats
#
# Default values:
#
RHYTHM ?= 4
BEATS  ?= 20

gen:
	@echo "Generating signal..."
	$(eval CMD := ./generate_signal)
	$(eval CMD += --rhythm $(RHYTHM))
	$(eval CMD += --beats $(BEATS))
	$(CMD)
	@echo "Done."

#
# Internal use:
#
libraries:
	g++ -std=c++20 -c cleaner.cpp
	g++ -std=c++20 -c visualize.cpp
