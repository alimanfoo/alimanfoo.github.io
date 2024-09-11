#!/bin/bash
while inotifywait -e close_write slides.ipynb
do
    pixi run jupyter nbconvert slides.ipynb \
	 --to slides \
	 --no-prompt \
	 --TagRemovePreprocessor.enabled=True \
	 --TagRemovePreprocessor.remove_input_tags remove_input \
	 --TagRemovePreprocessor.remove_all_outputs_tags remove_output \
	 --TagRemovePreprocessor.remove_cell_tags remove_cell \
         --stdout \
         > index.html
done
