"""This package provides classes and tools for the EMC3-EIRENE code.

This package provides classes that represents the EMC3-EIRENE input and 
output files for an easy access with python.
Furthermore, tools for analysing simulations obtained by EMC3-EIRENE are 
provided.

"""

from GeometryParameters import GeometryParameters
from PlasmaParameters import PlasmaParameters
from NeutralParameters import NeutralParameters
from ImpurityParameters import ImpurityParameters
from IterationControlParameters import IterationControlParameters
from Info import Info
from Zone import Zone
from Grid import Grid
from MagneticField import MagneticField
from FluxConservation import FluxConservation
from PhysicalCells import PhysicalCells
from PlasmaField import PlasmaField
from PlateCells import PlateCells
from Depo import Depo
from ExportToXdmf import ExportToXdmf
from TestCase import TestCase
from Installation import Installation
from AdditionalSurfaces import AdditionalSurfaces
from TargetFile import TargetFile
from TargetProfiles import TargetProfiles
from TargetProfiles import TargetProfilesElement
import Analysis
