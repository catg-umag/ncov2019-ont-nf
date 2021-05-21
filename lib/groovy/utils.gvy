import org.yaml.snakeyaml.Yaml

def addDefaultParamValues(params, defaultsFile) {
  yaml = new Yaml()
  defaultValues = yaml.load(file(defaultsFile))

  defaultValues.each { name, value ->
    params[name] = evaluateWorkflowVars(getParamValue(params, name, value))
  }
}

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
      switch(conf.validate) {
        case 'file':
          f = file(params[name])
          valid = f.exists() && f.isFile()
          break
        case 'directory':
          d = file(params[name])
          valid = d.exists() && d.isDirectory()
          break
        case 'boolean':
          valid = (params[name] instanceof Boolean)
          break
        case 'integer':
          valid = (params[name] instanceof Integer)
          break
      }

      if (!valid) {
        throwError("parameter ${name} must be a valid ${conf.validate}")
      }
    }
  }
}


def evaluateWorkflowVars(value) {
  matches = (value =~ /\$\{workflow\.[^}]*\}/).findAll()
  matches.each {
    name = (it =~ /(?<=\{workflow\.).*(?=\})/).findAll()[0]
    value = value.replace(it, workflow[name].toString())
  }

  return value
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
