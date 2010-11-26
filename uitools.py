import enthought.traits.api as eta
import enthought.traits.ui.api as etua
import simplejson as sj
import os

def strlist(x):
    return str(list(x))



################################################################################
class Workspace(eta.HasTraits):

    update = eta.Button("Fetch Variable")
    fname = eta.File
    v = ["banane"]
    ed = etua.SetEditor(values  = v, ordered = True, left_column_title  = 'to loose', right_column_title = 'to save')
    ordered_set = eta.List(editor = ed)
    tosave = eta.Button("Save")

    
    def __init__(self):
        """
        """

        self._update_who()


    def _update_who(self):
        """
        update value or content of %who
        """

        self.v = _ip.magic("who_ls")
        self.ed.values = self.v
        self.ordered_set.editor = self.ed

    
    def _update_fired(self):
        """
        update value or content of %who
        """


    def _tosave_fired(self):
        """
        write content of %who in choosen file
        """
        
        if self.fname == None:
            return 1
        else:
            print(self.ordered_set + "to save in " + self.fname)


    view = etua.View(
        etua.Item("update", show_label=False),
        etua.Item("ordered_set", style="custom", show_label=False),
        etua.Item("fname", style="simple", show_label=False),
        etua.Item("tosave", show_label=False),
        resizable = True
        )



def workspace():
    Workspace().configure_traits()


###################################################################################

class UiImport(eta.HasTraits):
    data = eta.Dict
    fname = eta.File("/Users/raoul/")
    openNanoQtf = eta.Button("Open NanoQt file")
    openTSVf = eta.Button("Open TSV file")
    
    def __init__(self, fname = ""):
        if fname != "":
            self.fname = fname

    def _openNanoQtf_fired(self):
        self.data = sj.load(file(self.fname, "r"))
        openNanoQtgraph(self.data)


    def _openTSVf_fired(self):
        openTSVgraph(self.fname)


    view = etua.View(
        etua.Item("data", style="simple"),
        etua.Item("fname", style="custom"),
        etua.Item("openNanoQtf", show_label = False),
        etua.Item("openTSVf", show_label = False),
        # resizable = True,
        scrollable = True,
        height = 640,
        width = 800
        )


def uiimport(fname = ""):
    wxuiimport = UiImport(fname)
    wxuiimport.configure_traits()
