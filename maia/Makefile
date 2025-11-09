SHELL := /bin/sh

# find build directories
BUILDDIRS = $(subst ./,,$(shell find -L . -mindepth 2 -maxdepth 2 -type f -name 'Makefile' -exec dirname '{}' \; \
						| sed '/doc/d' \
						| awk '!seen[$$0]++' \
						| echo $$(tr '\012' ' ')))

# defining canned recipes

define link_verbose =
ln -sf ./$(1)/bin/maia maia_$(subst build_,,$(1));
endef

define link_simple =
ln -sf ./$(1)/bin/maia maia;
endef

define remove_links =
rm -f ./maia*;
endef

define echo_target =
echo "... $(1)";
endef

# defining targets

# dynamic target generation: finds all targets defined for each build type
TOPTARGETS = $(shell sed -n '/^.PHONY/p;' $(foreach __bdir__,$(BUILDDIRS),$(__bdir__)/Makefile) \
	| sed 's/.PHONY\s\{,\}:\s\{,\}//g;/default_target/d;/distclean/d;/help/d' \
	| awk '!seen[$$0]++' \
	| echo $$(tr '\012' ' '))

DOCTARGETS = doc doc-fast

RMTARGETS = distclean

SUBDIRTARGETS = $(foreach __bdir__,$(BUILDDIRS),$(patsubst %,$(__bdir__)/%,$(TOPTARGETS)))

LINKTARGETS = $(foreach __bdir__,$(BUILDDIRS),$(patsubst %,%/link,$(__bdir__)))

.DEFAULT_GOAL = all

default_target : all

# Target execution recipes

$(TOPTARGETS) : $(BUILDDIRS)

$(BUILDDIRS) :
	@$(MAKE) -s -C ./$@ $(filter-out $@,$(filter-out link,$(MAKECMDGOALS)))

$(SUBDIRTARGETS) :
	@$(MAKE) -s -C ./$(wordlist 1,2,$(subst /, , $@))$(findstring /fast,$@)

$(DOCTARGETS) :
	@$(MAKE) -s -C ./doc $@

$(RMTARGETS) :
	@echo "Making $@..."
	@$(call remove_links)
	@rm -rf $(BUILDDIRS)
	@rm -f config.log
	@echo "done!"

.PHONY : $(BUILDDIRS) $(TOPTARGETS) $(SUBDIRTARGETS) $(RMTARGETS) $(DOCTARGETS)

# special link targets

link : $(BUILDDIRS)
	@$(call remove_links)
	@echo "Linking all applicable maia executables..."
	@$(foreach __dir__,$(BUILDDIRS),$(call link_verbose,$(__dir__)))

$(LINKTARGETS) :
	@$(MAKE) -s -C ./$(subst /link,,$@)
	@$(call remove_links)
	@echo "Linking the maia executable found in $(subst /link,,$@)."
	@$(call link_simple,$(subst /link,,$@))

.PHONY : link $(LINKTARGETS)

#Help Target
help:
	@echo "The following are some of the valid top level targets for this Makefile:"
	@echo "The default target is set to '$(.DEFAULT_GOAL)'"
	@$(foreach tar,$(TOPTARGETS),$(call echo_target,$(tar)))
	@echo "... link (calls the default target set to '$(.DEFAULT_GOAL)' and links all maia executables in Solver/"
	@echo "... $(DOCTARGETS) (builds the documentation in the doc/ directory)"
	@echo "The following top level targets remove all links and build directories:"
	@$(foreach tar,$(RMTARGETS),$(call echo_target,$(tar)))
	@echo "Specific build directories are recognized if they contain a Makefile and can be targeted individually"
	@echo "or be combined with a valid top level target by specifying <DIR>/<TARGET>."
	@echo "The following build-directory specific targets have been created:"
	@$(foreach tar,$(SUBDIRTARGETS),$(call echo_target,$(tar)))
	@echo "The following build-directory specific targets create a symbolic link to the maia executable in Solver"
	@echo "... $(foreach __bdir__,$(BUILDDIRS),$(__bdir__)/link)"
.PHONY : help
