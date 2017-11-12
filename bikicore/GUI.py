"""Shared classes for GUIs
""" 

from traits.api import HasTraits, Instance, List
from traitsui.api import View, Menu, MenuBar, Action, Handler, NoButtons
from traitsui.item import Item
import bikipy.bikicore.components as bkcc
import bikipy.bikicore.model as bkcm