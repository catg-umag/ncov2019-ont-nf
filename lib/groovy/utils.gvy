import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import org.yaml.snakeyaml.Yaml


def addDefaultParamValues(params, defaultsFile) {
  yaml = new Yaml()
  defaultValues = yaml.load(file(defaultsFile))

  defaultValues.each { name, value ->
    params[name] = evaluateWorkflowVars(getParamValue(params, name, value))
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

def getDateStr() {
  now = LocalDateTime.now()
  formatter = DateTimeFormatter.ofPattern('yyyyMMdd')

  return now.format(formatter)
}
