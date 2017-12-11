# named colors
COLOR_RST = \e[0m
COLOR_ACT = \e[1;32m
COLOR_ARG = \e[1;35m

# build tools and flags
CC = clang++
LD = clang++

# debug flags -D DEBUGGING_ENABLED -g
CCFLAGS = -fno-operator-names -march=native -std=c++14 -pthread -Wall -Wextra -fopenmp -fsanitize=undefined -D DEBUGGING_ENABLED -O3
LDFLAGS = -lfftw3 -pthread -Wl,-rpath -Wl,/usr/lib/openmpi -Wl,--enable-new-dtags -L/usr/lib/openmpi -lmpi_cxx -lmpi -fopenmp -fsanitize=undefined
OBJECTS = $(BUILD)/src/fft/fft.o $(BUILD)/src/grid.o $(BUILD)/src/main.o $(BUILD)/src/solver.o
BUILD = build
BIN = $(BUILD)/build

# phonies
.PHONY: all clean test reconf rebuild
all: $(BIN)
clean:
	@echo -e "$(COLOR_ACT)removing $(COLOR_ARG)$(BUILD)$(COLOR_RST)..."
	rm -rf $(BUILD)/
test: all
	@echo -e "$(COLOR_ACT)running $(COLOR_ARG)$(BIN)$(COLOR_RST)..."
	$(BIN)
reconf:
	@echo -e "$(COLOR_ACT)reconfiguring$(COLOR_RST)..."
	./configure
rebuild: clean
	@$(MAKE) --no-print-directory all

# build rules
$(BUILD)/:
	@echo -e "$(COLOR_ACT)making directory $(COLOR_ARG)$(BUILD)/$(COLOR_RST)..."
	mkdir -p $(BUILD)/
$(BUILD)/src: | $(BUILD)/
	@echo -e "$(COLOR_ACT)making directory $(COLOR_ARG)$(BUILD)/src$(COLOR_RST)..."
	mkdir -p $(BUILD)/src
$(BUILD)/src/fft: | $(BUILD)/src
	@echo -e "$(COLOR_ACT)making directory $(COLOR_ARG)$(BUILD)/src/fft$(COLOR_RST)..."
	mkdir -p $(BUILD)/src/fft
$(BUILD)/src/fft/fft.o: src/fft/fft.cc src/fft/fft.hh src/grid.hh src/utils/constant.hh src/utils/tools.hh src/utils/type.hh | $(BUILD)/src/fft
	@echo -e "$(COLOR_ACT)compiling $(COLOR_ARG)src/fft/fft.cc$(COLOR_RST)..."
	$(CC) -c -o '$@' '$<' $(CCFLAGS)
$(BUILD)/src/grid.o: src/grid.cc src/grid.hh src/utils/constant.hh src/utils/type.hh | $(BUILD)/src
	@echo -e "$(COLOR_ACT)compiling $(COLOR_ARG)src/grid.cc$(COLOR_RST)..."
	$(CC) -c -o '$@' '$<' $(CCFLAGS)
$(BUILD)/src/main.o: src/main.cc src/grid.hh src/solver.hh src/utils/constant.hh src/utils/type.hh | $(BUILD)/src
	@echo -e "$(COLOR_ACT)compiling $(COLOR_ARG)src/main.cc$(COLOR_RST)..."
	$(CC) -c -o '$@' '$<' $(CCFLAGS)
$(BUILD)/src/solver.o: src/solver.cc src/fft/fft.hh src/grid.hh src/solver.hh src/utils/constant.hh src/utils/type.hh | $(BUILD)/src
	@echo -e "$(COLOR_ACT)compiling $(COLOR_ARG)src/solver.cc$(COLOR_RST)..."
	$(CC) -c -o '$@' '$<' $(CCFLAGS)
$(BIN): $(OBJECTS) | $(BUILD)/
	@echo -e "$(COLOR_ACT)loading $(COLOR_ARG)build$(COLOR_RST)..."
	$(LD) -o '$@' $(OBJECTS) $(LDFLAGS)

