# System imports
import numpy
import sys
import unittest

from PyDakota import options_list
from PyDakota.options_list import *

class OptionsListTestCase(unittest.TestCase):
    "TestCase class for Teuchos.ParameterList"

    def test_set(self):
        opts = OptionsList()
        opts.set("key1",2.)
        assert opts.get("key1")==2.
        assert opts.get("key2")==None

        opts.set("key1",2)
        assert opts.get("key1")==2

        opts.set("key1","a")
        assert opts.get("key1")=="a"

        opts.set("key1",{})
        assert opts.get("key1")=={}

    def test_len(self):
        opts = OptionsList()
        assert len(opts)==0

        opts = OptionsList({'key1':'a','key2':2})
        assert len(opts)==2

    def test_pydict_to_options_list(self):
        opts = OptionsList()
        opts.set("key1",{})
        assert len(opts)==1

    def test__eq__(self):
        opts1 = OptionsList()
        opts1.set("key1",{})

        opts2 = OptionsList()
        opts2.set("key1",{})

        assert opts1==opts2

        assert opts1==opts1

        opts2.set("key2","a")
        assert opts1 != opts2

    def test__get_item__(self):
        opts = OptionsList()
        opts.set("key1",2.)
        assert opts["key1"]==2.

    def test__contains__(self):
        opts = OptionsList()
        opts.set("key1",2.)
        assert "key1" in opts

    def test__set_item__(self):
        opts = OptionsList()
        opts["key1"]=2.
        opts["key1"]==2.

    def test__str__(self):
        #print opts
        #print opts.__repr__
        pass

    def test_options_list_typemap_in(self):
        """
        Test that a wrapped function that takes an OptionsList as input
        also accepts a PythonDictionary
        """
        opts = OptionsList()
        str_types = ['int','double','string','optionslist']
        items = [1,2.,'a',OptionsList()]
        names = ['key%s'%(i+1) for i in range(len(items))]
        for i,type_str in enumerate(str_types):
            item = items[i]
            name = names[i]
            set_entry = options_list.__dict__[type_str + "set_entry"]
            tmp_opts = set_entry({},name,item)
            opts = set_entry(opts,name,item)
            assert tmp_opts=={names[i]:items[i]}

        pydict = dict((names[i],items[i]) for i in range(len(items)))
        assert opts == pydict

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(OptionsListTestCase))

    # Run the test suite
    result = unittest.TextTestRunner(verbosity=1).run(suite)

    
