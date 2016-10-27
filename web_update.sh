#!/bin/bash

cd cviceni
pdflatex mkp_cvika.ltx
pdflatex mkp_cvika.ltx
cp -f mkp_cvika.pdf ../web/pdf

cd ..
cd prednaska
pdflatex mkp_prednaska.tex
pdflatex mkp_prednaska.tex
cp -f mkp_prednaska.pdf ../web/pdf
