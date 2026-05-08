"use client";

import { useEffect, useState } from "react";

export type FetchState<T> =
  | { status: "loading" }
  | { status: "ready"; data: T }
  | { status: "error"; error: string };

export function useApi<T>(fetcher: () => Promise<T>, deps: unknown[] = []): FetchState<T> {
  const [state, setState] = useState<FetchState<T>>({ status: "loading" });
  useEffect(() => {
    let live = true;
    setState({ status: "loading" });
    fetcher()
      .then((data) => {
        if (live) setState({ status: "ready", data });
      })
      .catch((e: Error) => {
        if (live) setState({ status: "error", error: e.message });
      });
    return () => {
      live = false;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, deps);
  return state;
}
