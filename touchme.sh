for f in core/**/*.py; do bname="${f%.*}"; if [ "${bname: -2}" != __ ]; then touch "test/${bname}_test.py"; fi; done
