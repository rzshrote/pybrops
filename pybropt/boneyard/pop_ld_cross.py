import numpy

def pop_ld_cross(pRec, mpAB, fpAB, mfreq, ffreq):
    return ((0.5 * (1 - pRec) * (mpAB + fpAB)) +
            (0.5 * pRec * ((mfreq[:,None] * ffreq) + (ffreq[:,None] * mfreq))))
