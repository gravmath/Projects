(defun sentence    () (append (noun-phrase) (verb-phrase)))
(defun noun-phrase () (append (Article)     (Noun)))
(defun verb-phrase () (append (Verb)        (noun-phrase)))
(defun Article     () (one-of '(the a this that which these those whose)))
(defun Noun        () (one-of '(man ball woman table gamer snail book death poster wall paint)))
(defun Verb        () (one-of '(posted reviewed attacked punched hit took saw liked)))
(defun one-of      (set)
                   "Pick one element of a set and return it in a list"
                   (list (random-element set)))
(defun random-element (choices)
                      "Choose an element from a list at random"
                      (elt choices (random (length choices))))