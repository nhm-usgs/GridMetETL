#!/usr/bin/env python

"""Tests for `gridmetetl` package."""


import pytest
import numpy as np
from netCDF4 import default_fillvals

from gridmetetl import etl
from gridmetetl.helper import getaverage, np_get_wval


def test_getaverage():
    assert getaverage(np.array([2., 2.]), np.array([1., 1.])) == 2.0
    assert getaverage(np.array([2., 2.]), np.array([0.5, 0.5])) == 2.0
    assert getaverage(np.array([4., 2.]), np.array([0.75, 0.25])) == 3.5
    # test that array with nan return nan
    assert np.isnan(getaverage(np.array([2.0, np.nan]), np.array([1., 1.])))
    assert np.isnan(getaverage(np.array([np.nan, np.nan]), np.array([1., 1.])))

def test_np_get_wval():
    assert np_get_wval(np.array([2., 2.]), np.array([1., 1.])) == 2.0
    assert np_get_wval(np.array([2., 2.]), np.array([0.5, 0.5])) == 2.0
    assert np_get_wval(np.array([4., 2.]), np.array([0.75, 0.25])) == 3.5
    # test returns masked value if at least one value is nan
    assert np_get_wval(np.array([2.0, np.nan]), np.array([1., 1.])) == 2.0
    assert np_get_wval(np.array([2.0, np.nan]), np.array([0.5, 0.5])) == 2.0
    assert np_get_wval(np.array([4.0, np.nan]), np.array([0.75, 0.5])) == 4.0
    assert np_get_wval(np.array([4.0, 4.0, np.nan]), np.array([0.15, 0.75, 0.35])) == 4.0
    # test returns default value when array all nans
    assert np_get_wval(np.array([np.nan, np.nan]), np.array([1., 1.])) == default_fillvals['f8']

