import datetime
import sys
import os


def print_header(logger=None):

    msg = """
        ***********************************************************************
                           Generation of Conformers (GeCoS)
                         -----------------------------------
                         
                                    Version 0.1
                         
                                  Dr. Javier Ramos
                          Macromolecular Physics Department
                    Instituto de Estructura de la Materia (IEM-CSIC)
                                   Madrid (Spain)
                                   
                GeCoS is an open-source python library to quickly generate
                conformers of small molecules or polymer segments 
                using RdKit and OpenBabel libraries. Once, the
                conformers are generated, QM optimizations and subsequent 
                clustering can be done.
                
                This software is distributed under the terms of the
                GNU General Public License v3.0 (GNU GPLv3). A copy of 
                the license (LICENSE.txt) is included with this distribution. 
                
        ***********************************************************************
                     
        """

    print(msg) if logger is None else logger.info(msg)


# =============================================================================
def print_header_analysis(logger=None):
    msg = """
        ***********************************************************************
                           Generation of Conformers (GeCoS)
                                       Analysis
                         -----------------------------------

                                    Version 0.1

                                  Dr. Javier Ramos
                          Macromolecular Physics Department
                    Instituto de Estructura de la Materia (IEM-CSIC)
                                   Madrid (Spain)

                GeCoS is an open-source python library to quickly generate
                conformers of small molecules or polymer segments 
                using RdKit and OpenBabel libraries. Once, the
                conformers are generated, QM optimizations and subsequent 
                clustering can be done.

                This software is distributed under the terms of the
                GNU General Public License v3.0 (GNU GPLv3). A copy of 
                the license (LICENSE.txt) is included with this distribution. 

        ***********************************************************************

        """

    print(msg) if logger is None else logger.info(msg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger is None else logger.info(m)

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    m += "\t\t\t         or\n"
    m += "\t\t\tgecos_analysis".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger is None else logger.info(m)