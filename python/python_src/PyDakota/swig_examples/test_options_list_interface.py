#from OptionsList import OptionsList
from options_list_interface import *
def test_options_list_typemap_in():
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
    test_options_list_typemap_in()
