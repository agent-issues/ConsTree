#!/bin/bash
# TreeAnnotator 2.7.7 launcher with a configurable heap (the stock script fixes -Xmx8g).
B=/nobackup/pjjg18/constree-validation/ref/beast2.7.7/beast
exec "$B/jre/bin/java" -Dlauncher.wait.for.exit=true -Xss256m -Xmx${TA_XMX:-8g} \
  -Djava.library.path="$B/lib" -Duser.language=en -cp "$B/lib/launcher.jar" \
  beast.pkgmgmt.launcher.TreeAnnotatorLauncher "$@"
