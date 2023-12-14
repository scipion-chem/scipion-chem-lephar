================================
LePhar scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

`LePhar <http://www.lephar.com/index.htm>`_ software plugin defining objects and protocols for docking.

Full documentation to this plugin can be found in the `official documentation page <https://scipion-chem.github.io/docs/plugins/lephar/index.html>`_.

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

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-lephar.git

2. **Switch to the desired branch** (master or devel):

Scipion-chem-lephar is constantly under development and including new features.
If you want a relatively older an more stable version, use master branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-lephar
            git checkout devel

3. **Install**:

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
