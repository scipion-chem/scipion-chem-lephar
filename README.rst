================================
LePhar scipion plugin
================================

LePhar software plugin defining objects and protocols for docking.
You can check their website in http://www.lephar.com/index.htm

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version (Currently not available)**

.. code-block:: 

      scipion3 installp -p scipion-chem-lephar
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-lephar.git

2. Install:

.. code-block::

            scipion3 installp -p path_to_scipion-chem-lephar --devel

=================
Included programs
=================

Currently, the following programs are included in the Scipion interface:
    - **LePro**: For protein structure preparation
    - **LeDock**: For protein-ligand docking

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_dev.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/bioinformatics_prod.svg
