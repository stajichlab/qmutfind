// Saved results: a small key-value store in the browser's IndexedDB (database "qmutfind", store "kv"). Every call opens the database
// and closes it again, so no connection stays open. A browser can refuse IndexedDB (a private window, blocked site data), so
// every function can reject, and the page works without it. `?nostore=1` makes the page behave as if IndexedDB were missing
// (the tests use it).
const DATABASE = 'qmutfind';
const STORE = 'kv';

export function storeAvailable() {
  try {
    return typeof indexedDB !== 'undefined' && indexedDB !== null && !new URLSearchParams(location.search).has('nostore');
  } catch (error) {
    return false;
  }
}

function open() {
  return new Promise((resolve, reject) => {
    if (!storeAvailable()) { reject(new Error('IndexedDB is not available in this browser')); return; }
    let request;
    try { request = indexedDB.open(DATABASE, 1); } catch (error) { reject(error); return; }
    request.onupgradeneeded = () => request.result.createObjectStore(STORE);
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error || new Error('the browser could not open its store'));
    request.onblocked = () => reject(new Error('the browser store is blocked by another tab'));
  });
}

/**
 * Run `work(store, abort)` in one transaction; resolves with its result when the transaction has completed. `abort(error)` cancels
 * every change of the transaction and rejects with `error`. An error thrown by `work` does the same, so a write is all or none.
 */
async function transaction(mode, work) {
  const db = await open();
  try {
    return await new Promise((resolve, reject) => {
      const tx = db.transaction(STORE, mode);
      let result;
      let reason = null;
      const abort = (error) => { reason = error; try { tx.abort(); } catch (e) { /* already finished */ } };
      tx.oncomplete = () => resolve(result);
      tx.onerror = () => reject(reason || tx.error || new Error('the browser store failed'));
      tx.onabort = () => reject(reason || tx.error || new Error('the browser store was stopped'));
      try { result = work(tx.objectStore(STORE), abort); } catch (error) { abort(error); }
    });
  } finally {
    db.close();
  }
}

/** The values of the keys, as an object ({} for a key that is not there). */
export async function readStore(keys) {
  const found = {};
  await transaction('readonly', (store) => {
    for (const key of keys) {
      const request = store.get(key);
      request.onsuccess = () => { if (request.result !== undefined) found[key] = request.result; };
    }
  });
  return found;
}

/**
 * Write `entries` ({ key: value }) in one transaction: all or none. A value of `null` deletes its key. With `{ expect }` the write is
 * made only if the `saved` time of the stored `header` is `expect` (null when no header is stored): otherwise nothing is written and the
 * Error has `code` "CHANGED", because another tab saved in the meantime.
 */
export async function writeStore(entries, options = {}) {
  await transaction('readwrite', (store, abort) => {
    const put = () => {
      try {
        for (const [key, value] of Object.entries(entries)) {
          if (value === null) store.delete(key); else store.put(value, key);
        }
      } catch (error) {
        abort(error);
      }
    };
    if (!('expect' in options)) { put(); return; }
    const request = store.get('header');
    request.onsuccess = () => {
      const current = request.result ? request.result.saved : null;
      if (current === options.expect) { put(); return; }
      const error = new Error('the saved results were changed by another tab');
      error.code = 'CHANGED';
      abort(error);
    };
  });
}

/** Delete everything the page saved. */
export async function clearStore() {
  await transaction('readwrite', (store) => { store.clear(); });
}

/** { usage, quota } of this site's browser storage in bytes, or null when the browser cannot say. */
export async function estimateStorage() {
  try {
    if (navigator.storage && navigator.storage.estimate) {
      const { usage, quota } = await navigator.storage.estimate();
      return { usage, quota };
    }
  } catch (error) { /* not available */ }
  return null;
}
