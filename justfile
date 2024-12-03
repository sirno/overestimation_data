
create_figures:
  mkdir -p figures
  for clock in strict_clock relaxed_clock; do \
    for rate in nmr hmr; do \
      for pop in pop_tiny pop_small; do \
        python notebooks/overestimations.py --selector ${clock}/${rate}/${pop} --load-store; \
      done; \
    done; \
  done;
  for clock in strict_clock relaxed_clock; do \
    python notebooks/tables.py --selector ${clock}; \
  done;

list_xml:
  find data -name "template.xml" | sort
