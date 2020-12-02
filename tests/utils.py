"""
Utility functions for tests
"""
import unittest
from functools import wraps
import warnings

class FutureProofTestRunner(unittest.TextTestRunner):
    """
    A test runner that raises errors on Deprecation warnings.
    """
    def run(self, *args, **kwargs):
        warnings.filterwarnings("error", category=DeprecationWarning)
        return super().run(*args, **kwargs)

def retry_test(max_retries):
    """ Retry the decorated test up to `max_retries` times if failure happens. """

    def decorator(test_item):
        @wraps(test_item)
        def new_test_item(self, *args, **kwargs):
            tries = 1
            while tries < max_retries:
                try:
                    r = test_item(self, *args, **kwargs)
                except self.failureException:
                    tries += 1
                else:
                    return r

            # Final retry.
            return test_item(self, *args, **kwargs)

        return new_test_item

    return decorator
