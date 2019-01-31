
Changelog
=========

0.6.1
-----

* Fixed an issue where cache directories on MacOS could not be used.
* `CODParser` will now try to download files from three Crystallography Open Database mirrors.
* `AtomicStructure` and subclasses now support "truthiness", i.e. they are considered `False` if empty, and `True` otherwise.
* Added the `AtomicStructure.satisfying` method, to extract atoms satisfying a predicate from structures
* Added the `is_element` function. It can be used to make `AtomicStructure.satisfying` more readable.