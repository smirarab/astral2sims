#!/bin/bash

grep finish model.*/*/*.astral.err | sed -e "s/\// /g" -e "s/model.//g" -e "s/\./ /" -e "s/000\./K /" -e "s/000K /M /" -e "s/\.genes/ /" -e "s/.astral.err:.*in / /g"|sed -e "s/-halfresolved//g"| tee runtime.stat
cat mrlruntime | sed -e "s/\// /g" -e "s/model.//g" -e "s/\./ /" -e "s/000\./K /" -e "s/000K /M /" -e "s/\.genes/ /" -e "s/.halfresolved//g" -e "s/$/ secs/" | tee -a runtime.stat
cat greedytime | sed -e "s/\// /g" -e "s/model.//g" -e "s/\./ /" -e "s/000\./K /" -e "s/000K /M /" -e "s/\.genes/ /" -e "s/.halfresolved//g" -e "s/$/ secs/" | tee -a runtime.stat
grep real model.*/*/*.time.stat| sed -e "s/\// /g" -e "s/model.//g" -e "s/\./ /" -e "s/000\./K /" -e "s/000K /M /" -e "s/\.genes/ /" -e "s/.time.stat:real / /g"|sed -e "s/.halfresolved//g"|tee -a runtime.stat
grep Overall model.*/*/RAxML_info.mrl.*|sed -e "s/\// /g" -e "s/model.//g" -e "s/\./ /" -e "s/000\./K /" -e "s/000K /M /" -e "s/\.genes/ /" -e "s/ RAx.*mrl/ mrl/g" -e "s/:Over.*: / /" -e "s/ or.*//"|sed -e "s/.halfresolved//g" |tee -a runtime.stat
