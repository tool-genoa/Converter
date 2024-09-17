"""
    This file contains utility functions that are used in the converter
"""

def isfloat(value):
    """Check if a string is a float"""
    try:
        float(value)
        return True
    except ValueError:
        return False
