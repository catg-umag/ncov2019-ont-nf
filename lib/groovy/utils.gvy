import org.yaml.snakeyaml.Yaml

/*
 * Add default parameter values (a mapping) if they are not defined,
 * reading them from a YAML file
 */
def addDefaultParamValues(params, defaultsFile) {
  yaml = new Yaml()
  defaultValues = yaml.load(file(defaultsFile))

  defaultValues.each { name, value ->
    params[name] = evaluateWorkflowVars(getParamValue(params, name, value))
  }
}

/*
 * Check if the parameter values are valid (mandatory, conditional and type),
 * reading the validation rules from a YAML file
 */
def validateParameters(params, paramConfigFile) {
  yaml = new Yaml()
  paramConfig = yaml.load(file(paramConfigFile))

  paramConfig.each { name, conf ->
    // check required
    if (conf.required && isEmpty(params, name)) {
      throwError("parameter ${name} is required")
    }

    // check required_if
    if (
      conf.containsKey('required_if')
      && !isEmpty(params, conf.required_if)
      && isEmpty(params, name)
      ) {
      throwError("parameter ${name} is required if ${conf.required_if} is specified")
    }

    // check type
    if (!isEmpty(params, name) && conf.containsKey('validate')) {
      valid = true
      switch (conf.validate) {
        case 'file':
          f = file(params[name])
          valid = f.exists() && f.isFile()
          break
        case 'directory':
          d = file(params[name])
          valid = d.exists() && d.isDirectory()
          break
        case 'boolean':
          valid = (params[name].getClass() == Boolean)
          break
        case 'integer':
          valid = (params[name].getClass() == Integer)
          break
      }

      if (!valid) {
        throwError("parameter ${name} must be a valid ${conf.validate}")
      }
    }
  }
}

/*
 * Evaluate workflow instrospection variables (workflow.*) inside a string
 */
def evaluateWorkflowVars(value) {
  matches = (value =~ /\$\{workflow\.[^}]*\}/).findAll()
  evaluated = value
  matches.each { m ->
    name = (m =~ /(?<=\{workflow\.).*(?=\})/).findAll()[0]
    evaluated = evaluated.replace(m, workflow[name].toString())
  }

  return evaluated
}

/*
 * Check if file exists of throw error
 */
def pathCheck(path, isDirectory = false) {
  f = file(path)
  if (!f.exists()) {
    throwError("Path ${path} does not exist")
  } else if (!f.isFile() && !isDirectory) {
    throwError("Path ${path} is not a file")
  } else if (!f.isDirectory() && isDirectory) {
    throwError("Path ${path} is not a directory")
  }

  return f
}

def isEmpty(params, name) {
  return !params.containsKey(name) || params[name] == null
}

def getParamValue(params, name, defaultValue) {
  return isEmpty(params, name) ? defaultValue : params[name]
}

def throwError(msg) {
  exit 1, "ERROR: $msg"
}
