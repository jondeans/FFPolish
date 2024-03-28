VENVDIR=.venv
VENV=$(VENVDIR)/bin
VENV_PYTHON=$(VENV)/python
VENV_PIP=$(VENV)/pip

venv: $(VENV)/activate
$(VENV)/activate: pyproject.toml
	test -d $(VENVDIR) || python -m venv $(VENVDIR)
	$(VENV_PIP) install -U pip
	touch $(VENV)/activate

.PHONY: install
install: venv
	$(VENV_PIP) install -e .

.PHONY: install-dev
install-dev: venv
	$(VENV_PIP) install -e ".[dev]"

.PHONY: check
check: install-dev
	ruff check src/ tests/
	ruff format --check src/ tests/

.PHONY: clean
clean: venv
	rm -rf $(VENV_DIR)
	rm -rf *.egg-info

.PHONY: fix
fix: venv
	ruff src/ tests/
	ruff format src/ tests/

.PHONY: test
test: venv
	$(PYTHON) -m pytest tests
