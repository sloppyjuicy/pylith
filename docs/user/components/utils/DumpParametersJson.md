# DumpParametersJson

% WARNING: Do not edit; this is a generated file!
:Full name: `pylith.utils.DumpParametersJson`
:Journal name: `dumpparamters`

Dump PyLith parameter information to an ASCII file.

If you do not set the filename for the progress monitor, then PyLith will create one using the
simulation name from the application defaults settings.

Implements `DumpParameters`.

## Pyre Properties

* `filename`=\<str\>: Name of file written with parameters.
  - **default value**: ''
  - **current value**: '', from {default}
* `indent`=\<int\>: Nmber of spaces to indent, use a negative number for no newlines.
  - **default value**: 4
  - **current value**: 4, from {default}
* `style`=\<str\>: Style of JSON file [compact, normal].
  - **default value**: 'normal'
  - **current value**: 'normal', from {default}
  - **validator**: (in ['normal', 'compact'])

## Example

Example of setting `DumpParametersJson` Pyre properties and facilities in a parameter file.

:::{code-block} cfg
[pylithapp]
dump_parameters = pylith.utils.DumpParametersJson

[pylithapp.dump_parameters]
filename = output/parameters.json
style = normal
verbose = True
:::

