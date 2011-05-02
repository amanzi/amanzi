;;
;; cpplint.el
;;
;; Support for the Google code checking tool cpplint
;;
;; To load:
;;
;; Load this file into emacs. For example, add the following to your .emacs or .emacs.d/init.el:
;;
;;   (load-file "path/to/this/file/cpplint.el")
;;
;; You will need cpplint.py somewhere in your path.
;;
;; To activate in a particular buffer:
;;
;;   M-x activate-cpplint
;;
;; To deactivate in a particular buffer:
;;
;;  M-x deactivate-cpplint
;;
;; To use:
;;
;;   Warnings returned by cpplint will be highlighted. Use M-p, M-n to
;;   move between warnings and see the text of the warning in the
;;   message window.
;;
;;   Hover over a highlighted error with the mouse to see the message
;;   in a tooltip.
;;
;; Please report errors to mwbuksas@lanl.gov
;;

(require 'google-c-style)

(eval-after-load "flymake"
  '(progn (add-to-list 'flymake-allowed-file-name-masks
                       '("\\.hh\\'" flymake-cpplint-init flymake-simple-cleanup flymake-get-real-file-name))
          (add-to-list 'flymake-allowed-file-name-masks
                       '("\\.cc\\'" flymake-cpplint-init flymake-simple-cleanup flymake-get-real-file-name))
          (add-to-list 'flymake-allowed-file-name-masks
                       '("\\.cpp\\'" flymake-cpplint-init flymake-simple-cleanup flymake-get-real-file-name))
          (add-to-list 'flymake-allowed-file-name-masks
                       '("\\.hpp\\'" flymake-cpplint-init flymake-simple-cleanup flymake-get-real-file-name))))

(require 'flymake)


(defun flymake-cpplint-init ()
  (let* ((temp-file (flymake-init-create-temp-buffer-copy
                     'flymake-create-temp-inplace))
         (local-file (file-relative-name
                      temp-file
                      (file-name-directory buffer-file-name))))
    (list "cpplint.py" (list "--filter=-whitespace,-legal/copyright" local-file))))

;; From http://www.emacswiki.org/emacs/FlyMake

(defvar my-flymake-minor-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map "\M-p" 'flymake-goto-prev-error)
    (define-key map "\M-n" 'flymake-goto-next-error)
    map)
  "Keymap for my flymake minor mode.")

(defun my-flymake-err-at (pos)
  (let ((overlays (overlays-at pos)))
    (remove nil
            (mapcar (lambda (overlay)
                      (and (overlay-get overlay 'flymake-overlay)
                           (overlay-get overlay 'help-echo)))
                    overlays))))

(defun my-flymake-err-echo ()
  (message "%s" (mapconcat 'identity (my-flymake-err-at (point)) "\n")))

(defadvice flymake-goto-next-error (after display-message activate compile)
  (my-flymake-err-echo))

(defadvice flymake-goto-prev-error (after display-message activate compile)
  (my-flymake-err-echo))

(define-minor-mode my-flymake-minor-mode
  "Simple minor mode which adds some key bindings for moving to the next and previous errors.

Key bindings:

\\{my-flymake-minor-mode-map}"
  nil
  nil
  :keymap my-flymake-minor-mode-map)


(defun activate-cpplint () "Activate cpplint on the current buffer" (interactive)
  (progn (flymake-mode t)
         (my-flymake-minor-mode t)))

(defun deactivate-cpplint () "Deactivate cppline on the current buffer" (interactive)
  (progn (flymake-mode nil)
         (my-flymake-minor-mode nil)))


(provide 'cpplint)
