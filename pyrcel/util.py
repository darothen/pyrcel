"""
Software utilities

"""


class ParcelModelError(Exception):
    """Custom exception to throw during parcel model execution."""

    def __init__(self, error_str):
        self.error_str = error_str

    def __str__(self):
        return repr(self.error_str)
