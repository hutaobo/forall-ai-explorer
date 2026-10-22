// Delegated listener survives Streamlit rerenders without accumulating handlers.
(() => {
  const doc = window.parent.document;
  if (doc.forallPromptCompletion) doc.removeEventListener('keydown', doc.forallPromptCompletion, true);
  const complete = (event) => {
    const field = event.target;
    if (event.key !== 'Tab' || event.shiftKey || event.ctrlKey || event.altKey ||
        event.metaKey || event.isComposing || field.tagName !== 'TEXTAREA' ||
        field.getAttribute('aria-label') !== 'Ask AI') return;
    const example = field.placeholder;
    const typed = field.value;
    if (!example || typed === example ||
        !example.toLowerCase().startsWith(typed.toLowerCase()) ||
        field.selectionStart !== typed.length || field.selectionEnd !== typed.length) return;
    event.preventDefault();
    event.stopImmediatePropagation();
    const view = field.ownerDocument.defaultView;
    Object.getOwnPropertyDescriptor(view.HTMLTextAreaElement.prototype, 'value').set.call(field, example);
    field.dispatchEvent(new view.Event('input', {bubbles: true}));
    field.dispatchEvent(new view.Event('change', {bubbles: true}));
    field.setSelectionRange(example.length, example.length);
  };
  doc.forallPromptCompletion = complete;
  doc.addEventListener('keydown', complete, true);
})();
